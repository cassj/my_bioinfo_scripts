require 'rubygems'
require 'fileutils'
require 'AWS'
require 'json'


unless Capistrano::Configuration.respond_to?(:instance)
  abort "capistrano/ec2manager requires Capistrano 2"
end

Capistrano::Configuration.instance(:must_exist).load do
  
  namespace 'EC2' do

    def _cset(name, *args, &block)
      unless exists?(name)
        set(name, *args, &block)
      end
    end
    
    configuration.load do

      # User AWS details
      _cset (:access_key) {abort "Please specify your amazon access key, set :access_key, 'username'"}
      _cset (:secret_access_key) {abort "Please specify your amazon secret access key, set :secret_access_key, 'username'"}

      #default to EU West.
      _cset :ec2_url,  'eu-west-1.ec2.amazonaws.com'
      
      # localhost as master node by default
      _cset :master, 'localhost'
      
      #And create ec2 object
      set :ec2, AWS::EC2::Base.new( :access_key_id => ACCESS_KEY_ID, :secret_access_key => SECRET_ACCESS_KEY, :server => server)

    end
    
    
    ###################################################
    # Methods
    
    def load_instances(ridfile, ec2)
      
      unless ridfile
        puts "No ridfile specified in load_instances"
        exit 1
      end
      
      unless File.exists?(ridfile)
        puts "File #{ridfile} not found"
        exit 1
      end
      
      ridfile = File.open(ridfile, "r")
      unless ridfile
        puts "Can't open file #{ridfile} for reading"
        exit 1
      end
      
      rid = ridfile.gets.chomp
      ridfile.close
      
      inst  = ec2.describe_instances.reservationSet.item.select {|i| i.reservationId == rid }
      
      if inst.length == 0
        raise "No instances associated with reservation ID found. Please manually check your EC2 account."
      elsif inst.length >1
        puts "Multiple reservation sets associated with this reservation ID. Something has gone awry. Please manually check your EC2 account."
      end
      
      return inst[0].instancesSet.item
      
    end
    
    def instance_ids( insts )
      
      return insts.map{ |inst| inst.instanceId }
      
    end
    
    def instance_status( insts )
    
      stats  = Hash[]
      insts.each{ |inst| stats[inst.instanceId]  = inst.instanceState  }
      return stats
      
    end
    
    def instance_uri( insts )
      
      stats  = Hash[]
      insts.each{ |inst| stats[inst.instanceId]  = inst.dnsName }
      return stats
      
    end
    
    def instance_ami( insts )
      
      return insts.map{ |inst| inst.imageId }
      
    end



    ######################################
    # Tasks
    
    desc 'Start a group of EC2 images'
    task :start, :hosts => master do
      
      #user-supplied settings
      ridfile       = variables[:ridfile] || {abort "Specify an rid (reservation id) file to keep track of the instances you're using. For example cap start_ec2 -S ridfile='<filename>'"}
      nhosts        = variables[:nhosts] || 1
      ami           = variables[:ami] || {abort "Please specify an ami to use: cap start_ec2 -S ami='ami123456' "}
      instance_type = variables[:instance_type] || {abort "Please specify an instance_type to use: cap start_ec2 -S instance_type='m1.small'"}
      key           = variables[:key] || {abort "Please specify an AWS key to use: cap start_ec2 -S key='mykey' "}
      
      if File.exists?(ridfile)
        puts "File #{ridfile} already exists. Trying to use existing instances..."
        insts = load_instances(ridfile, ec2)
        
        #check number 
        unless insts.length == nhosts
          abort "Requested number of instances does not match running instances for this ridfile. Giving up."
        end
        
        #check ami
        amis = instance_ami( insts )
        unless amis.all?{|a| a == ami}
          abort "Instances for this ridfile do not match the requested AMI. Giving up "
        end
        
        #check status
        stats = instance_status( insts )
        codes = stats.keys.map{|k| stats[k].code}
        unless codes.all?{|c| Integer(c) == 16}
          abort "Some instances for this ridfile are not currently running. Giving up. "
        end
        
        puts "Success. Using existing instances from ridfile #{ridfile}"
        exit 0
        
      else
        ridfileFile = File.new(ridfile, "w")
      end
      
      unless ridfileFile
        abort "Can't create file #{ridfile}. for writing "
      end
      
      #type_info  = JSON.generate(instance_type_info[instance_type])
      
      #start required instances.
      new_inst = ec2.run_instances(
                                   :image_id => ami, 
                                   :min_count=> nhosts, 
                                   :max_count => nhosts, 
                                   :key_name => key,
                                   :instance_type => instance_type
                                   )
      #write the reservation Id out to a file so we can keep track of the machines we're using for these tasks.
      ridfileFile.syswrite(new_inst.reservationId)
      ridfileFile.close
      
      #block until they're up and running (or broken) so other tasks can't send stuff to them before they're ready
      ok_go = false
      
      while(!ok_go)
        
        puts "Instances pending, please wait ..."
        sleep(10)
        
        insts = load_instances(ridfile, ec2)
        stats = instance_status(insts)
        codes = stats.keys.map{|k| stats[k].code}
        
        #0 is pending, 16 is running, anything bigger is shutting down.
        if codes.any?{|c| Integer(c) > 16} 
          puts "Some instances seem to have stopped. Please check manually. \nShould be safe to run ec2_stop -S ridfile=#{ridfile} to shutdown remaining instances"
          exit 1
        end
        
        ok_go = codes.all?{ |c| Integer(c) == 16 }
        
      end
      
      puts "Started #{nhosts} instances of AMI #{ami}. ridfile #{ridfile}"
      
    end
    
    
    desc 'Stop the group of EC2 images'
    task :stop, :hosts => master do
      
      ridfile = variables[:ridfile]
      puts ridfile
      if ridfile.nil?
        puts "Specify an rid (reservation id) file correponding to the reservation set of instances you want to stop. For example cap stop_ec2 -S ridfile='<filename>'"
        exit 1
      end
      
      insts = load_instances(ridfile, ec2)
      ids = instance_ids(insts)
      
      ec2.terminate_instances(:instance_id => ids)
      
      sleep(10)
      insts = load_instances(ridfile, ec2)
      stats = instance_status(insts)
      codes = stats.keys.map{|k| stats[k].code}
      if codes.all?{|c| Integer(c) > 16}
        puts "Stopped instances with ids "
        puts stats.keys
        FileUtils.rm( ridfile ) 
      else
        puts "Failed to stop instances, please check manually"
        puts "Details:"
        puts "\t#{stats}"
        exit 1
      end
    end
    
    
    def claim( host )
      
      #return insts.map{ |inst| inst.imageId }
      
    end
    
    
  end

end



